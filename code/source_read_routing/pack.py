"""Lossless differential storage of full executions, with no event elision.

Later layers store the uint64 differences from the first layer's explicit
rows. Intervention runs store differences from the corresponding baseline
rows. All differences are retained, including metadata differences. Addition
modulo 2^64 recovers every original field, without consulting the simulator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from tape import decode, encode, hash_file, load_json


class Reader:
    def __init__(self, blocks):
        self.blocks = iter(blocks)
        self.pending = b""

    def read(self, size):
        pieces = []
        while size:
            if not self.pending:
                self.pending = next(self.blocks, b"")
                if not self.pending:
                    break
            part, self.pending = self.pending[:size], self.pending[size:]
            pieces.append(part)
            size -= len(part)
        return b"".join(pieces)


class Segment:
    def __init__(self, source, events, base=None, copy=None):
        self.source, self.remaining, self.base, self.copy = source, events*64, base, copy
        self.digest = hashlib.sha256()

    def read(self, size):
        if not self.remaining:
            return b""
        data = self.source.read(min(size, self.remaining))
        if not data or len(data)%64:
            raise ValueError("truncated event segment")
        self.remaining -= len(data)
        self.digest.update(data)
        if self.copy is not None:
            self.copy.write(data)
        if self.base is None:
            return data
        reference = self.base.read(len(data))
        if len(reference) != len(data):
            raise ValueError("delta reference length differs")
        return (np.frombuffer(data, dtype="<u8") - np.frombuffer(reference, dtype="<u8")).tobytes()


def segment_rows(folder, manifest_name, active=()):
    if manifest_name in active or Path(manifest_name).name != manifest_name or any(c in manifest_name for c in "\\/:"):
        raise ValueError("cyclic/nonlocal segment reference")
    manifest = load_json(folder/manifest_name)
    if manifest["schema"] != "oph.source_read_routing.segment.v1":
        raise ValueError("segment schema")
    if manifest.get("delta_base") in active+(manifest_name,):
        raise ValueError("cyclic segment reference")
    # Parts are fixed byte slices of one compressed stream, not event cuts.
    def parts():
        for part in manifest["parts"]:
            if Path(part["path"]).name != part["path"] or any(c in part["path"] for c in "\\/:"):
                raise ValueError("nonlocal tape part")
            path = folder/part["path"]
            if path.stat().st_size != part["bytes"] or hash_file(path) != part["sha256"]:
                raise ValueError("tape part custody")
            with path.open("rb") as stream:
                while data := stream.read(4*1024*1024):
                    yield data
    baseline = Reader(segment_rows(folder, manifest["delta_base"], active+(manifest_name,))) if manifest.get("delta_base") else None
    digest, size = hashlib.sha256(), 0
    for data in decode(Reader(parts())):
        if baseline is not None:
            reference = baseline.read(len(data))
            if len(reference) != len(data):
                raise ValueError("delta base length")
            data = (np.frombuffer(data, dtype="<u8")+np.frombuffer(reference, dtype="<u8")).tobytes()
        size += len(data); digest.update(data)
        yield data
    if baseline is not None and baseline.read(1):
        raise ValueError("delta base suffix")
    if size != manifest["events"]*64 or digest.hexdigest() != manifest["decoded_sha256"]:
        raise ValueError("segment decoded custody")


def pack(receipt_path, destination, workspace):
    original = json.loads(receipt_path.read_text(encoding="utf-8"))
    q, variant = original["inputs"]["q"], original["variant"]
    n, c, k = original["inputs"]["sites"], original["inputs"]["carriers"], original["inputs"]["rounds"]
    prelude = 13*c + 8*n
    layer_events, remainder = divmod(original["costs"]["events"]-prelude, k)
    if remainder:
        raise ValueError("inconsistent layer sizes")
    destination.mkdir(parents=True, exist_ok=True)
    workspace.mkdir(parents=True, exist_ok=True)
    source_input = receipt_path.parent/f"q{q}.input"
    shutil.copyfile(source_input, destination/source_input.name)
    names = []
    first_layer = workspace/f"q{q}.first_layer.raw"
    with (receipt_path.parent/original["tape"]["path"]).open("rb") as stream:
        reader = Reader(decode(stream))
        for phase, events in enumerate([prelude]+[layer_events]*k):
            name = f"q{q}_{variant}_phase{phase}"
            delta_name = None
            base, copy = None, None
            if variant != "baseline":
                delta_name = f"q{q}_baseline_phase{phase}.segment.json"
                base = Reader(segment_rows(destination, delta_name))
            elif phase > 1:
                delta_name = f"q{q}_baseline_phase1.segment.json"
                base = first_layer.open("rb")
            elif phase == 1:
                copy = first_layer.open("wb")
            segment = Segment(reader, events, base, copy)
            encoded_path = workspace/f"{name}.packed"
            with encoded_path.open("wb") as output:
                encoded = encode(segment, output)
            if hasattr(base, "close"):
                base.close()
            if copy is not None:
                copy.close()
            if encoded["events"] != events:
                raise ValueError("packing lost events")
            manifest = {"schema": "oph.source_read_routing.segment.v1", "events": events,
                        "decoded_sha256": segment.digest.hexdigest(), "parts": []}
            if delta_name:
                manifest["delta_base"] = delta_name
            with encoded_path.open("rb") as encoded_stream:
                part = 0
                while data := encoded_stream.read(32*1024*1024):
                    part_path = destination/f"{name}.{part:03d}.tape"
                    part_path.write_bytes(data)
                    manifest["parts"].append({"path": part_path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
                    part += 1
            manifest_path = destination/f"{name}.segment.json"
            manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2)+"\n", encoding="utf-8", newline="\n")
            names.append({"path": manifest_path.name, "sha256": hash_file(manifest_path)})
            print(f"Packed {name}: {events} explicit events, {sum(p['bytes'] for p in manifest['parts'])} bytes", flush=True)
        if reader.read(1):
            raise ValueError("unpacked event suffix")
    result = dict(original)
    result["tape"] = {"events": original["tape"]["events"], "decoded_bytes": original["tape"]["decoded_bytes"],
                      "decoded_sha256": original["tape"]["decoded_sha256"], "segments": names}
    out = destination/receipt_path.name
    out.write_text(json.dumps(result, sort_keys=True, indent=2)+"\n", encoding="utf-8", newline="\n")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    pack(args.receipt, args.destination, args.work)
