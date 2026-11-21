"""Lossless comparison of the prepared detector arrays and reported readings."""
import hashlib
import json
from pathlib import Path

import numpy as np

from tools.strict_json import load
from .kernel import require

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


# Digest of the input registration at 50c6d4fd9.
FROZEN_MANIFEST = "aa48fab83b02bd5b6047e2c02f734fd6ac8fcb377134ff1e773ea76cd97e1017"


def verify_freeze():
    manifest = HERE / "freeze.json"
    digest = hashlib.sha256(manifest.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(digest == FROZEN_MANIFEST, "changed frozen registration")
    freeze = load(manifest)
    for directory, pins in ((ROOT, freeze["simulator_files"]), (HERE, freeze["reviewed_fixtures_sha256"])):
        for name, expected in pins.items():
            path = directory / name
            raw = path.read_bytes()
            if path.suffix != ".npz":
                raw = raw.replace(b"\r\n", b"\n")
            require(hashlib.sha256(raw).hexdigest() == expected, f"changed frozen input: {name}")


def serialize(value):
    if isinstance(value,dict):
        return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,list):
        return [serialize(v) for v in value]
    if isinstance(value,(float,np.floating)):
        return float(format(value,'.9g'))
    return value


def receipt(path,value):
    content = json.dumps(serialize(value),indent=1,sort_keys=True,allow_nan=False)+'\n'
    stored = json.dumps(load(path),indent=1,sort_keys=True,allow_nan=False)+'\n'
    require(stored==content,f'changed result: {path.name}')


def archive(path,arrays):
    with np.load(path,allow_pickle=False) as stored:
        require(sorted(stored.files)==sorted(arrays),'changed detector menu')
        for k,a in arrays.items():
            require(stored[k].dtype==a.dtype and np.array_equal(stored[k],a),f'changed detector array: {k}')
        return dict(stored)
