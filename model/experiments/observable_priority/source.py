"""Read canonical transcripts and verify the source inputs for the experiment."""

from functools import cache
import gzip
import hashlib
import json
from pathlib import Path
import sys

from model.executable_model import carrier

from .observed import need


PIN = "6c3c142d3280b9a90f3b072fa5cc0be7581d226e"
EXPECTED = Path(__file__).with_name("expected.json")
ROOT = Path(carrier.__file__).resolve().parents[2]


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=True, separators=(",", ":"), allow_nan=False)


def fingerprint(value):
    return hashlib.sha256(canonical(value).encode("ascii")).hexdigest()


def read(path):
    def unique(pairs):
        obj = {}
        for key, value in pairs:
            need(key not in obj, "duplicate JSON field")
            obj[key] = value
        return obj

    def nonfinite(_):
        raise ValueError("nonfinite JSON number")

    path = Path(path)
    data = path.read_bytes()
    if path.suffix == ".gz":
        data = gzip.decompress(data)
    return json.loads(data, object_pairs_hook=unique, parse_constant=nonfinite)


@cache
def bind_source():
    expected = read(EXPECTED)
    need(expected["replay"]["source_commit"] == PIN, "source revision mismatch")
    observed = {}
    reference = expected["source_sha256"]
    for name in reference:
        path = ROOT / name
        data = path.read_bytes()
        if path.suffix == ".py":
            # Git can check out text files with CRLF on Windows.
            data = data.replace(b"\r\n", b"\n")
            module = ".".join(name[:-3].removesuffix("/__init__").split("/"))
            loaded = sys.modules.get(module)
            if loaded is not None:
                need(Path(loaded.__file__).resolve() == path.resolve(), "foreign module: " + module)
        observed[name] = hashlib.sha256(data).hexdigest()
    profiles = {PIN: reference}
    profiles.update({commit: reference | changes for commit, changes in expected["source_overrides"].items()})
    for commit, profile in profiles.items():
        if observed == profile:
            return expected, commit
    changed = [name for name in reference if observed[name] != reference[name]]
    raise ValueError("unrecognized source inputs: " + ", ".join(changed))


def board_hash(level):
    expected, _ = bind_source()
    return expected["source_sha256"][f"model/executable_model/carriers/FEDERATION_LEVEL{level}.npz"]
