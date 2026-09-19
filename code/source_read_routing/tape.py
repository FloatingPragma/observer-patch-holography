"""Lossless byte-plane compression of explicit 64-byte primitive-event rows.

No event, field, or value is inferred from a routing recipe. Each decoded row
contains all eight uint64 fields written by the executing native producer.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import struct
import zlib

import numpy as np

MAGIC = b"OPHW12T2"
ROWS = 65536
WIDTH = 64


def load_json(path):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result
    def reject(value):
        raise ValueError("nonintegral JSON number: " + value)
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique,
                      parse_float=reject, parse_constant=reject)


def encode(source, destination):
    digest = hashlib.sha256()
    count = 0
    destination.write(MAGIC)
    while True:
        chunks = []
        size = 0
        while size < ROWS * WIDTH:
            part = source.read(ROWS * WIDTH - size)
            if not part:
                break
            chunks.append(part)
            size += len(part)
        if not size:
            break
        data = b"".join(chunks)
        if size % WIDTH:
            raise ValueError("partial primitive event")
        digest.update(data)
        n = size // WIDTH
        rows = np.frombuffer(data, dtype="<u8").reshape(n, 8)
        ops = rows[:, 0].astype(np.uint8)
        if np.any(rows[:, 0] > 255):
            raise ValueError("opcode exceeds codec byte")
        permutation = np.argsort(ops, kind="stable")
        ordered = rows[permutation]
        differences = ordered.copy()
        differences[1:] -= ordered[:-1]
        zigzag = (differences << np.uint64(1)) ^ (differences.view("<i8") >> 63).view("<u8")
        planes = ops.tobytes() + zigzag.view(np.uint8).reshape(n, WIDTH).T.copy().tobytes()
        compressed = zlib.compress(planes, level=1)
        destination.write(struct.pack("<II", n, len(compressed)))
        destination.write(compressed)
        count += n
    return {"events": count, "decoded_sha256": digest.hexdigest(), "decoded_bytes": count * WIDTH}


def decode(source):
    magic = source.read(8)
    if magic not in (MAGIC, b"OPHW12T1"):
        raise ValueError("event tape magic")
    partial_seen = False
    while header := source.read(8):
        if len(header) != 8 or partial_seen:
            raise ValueError("truncated or noncanonical tape block")
        n, size = struct.unpack("<II", header)
        if not 0 < n <= ROWS or not 0 < size <= ROWS * WIDTH + 65536:
            raise ValueError("tape block bounds")
        partial_seen = n < ROWS
        compressed = source.read(size)
        if len(compressed) != size:
            raise ValueError("truncated tape payload")
        decoder = zlib.decompressobj()
        expected = n * (WIDTH + (magic == MAGIC))
        planes = decoder.decompress(compressed, expected + 1)
        if (len(planes) != expected or not decoder.eof or decoder.unused_data
                or decoder.unconsumed_tail):
            raise ValueError("invalid tape compression")
        if magic == b"OPHW12T1":
            yield np.frombuffer(planes, dtype=np.uint8).reshape(WIDTH, n).T.copy().tobytes()
        else:
            ops = np.frombuffer(planes[:n], dtype=np.uint8)
            zigzag = np.frombuffer(planes[n:], dtype=np.uint8).reshape(WIDTH, n).T.copy().view("<u8").reshape(n, 8)
            differences = (zigzag >> np.uint64(1)) ^ (np.uint64(0) - (zigzag & np.uint64(1)))
            ordered = np.cumsum(differences, axis=0, dtype=np.uint64)
            rows = np.empty_like(ordered)
            rows[np.argsort(ops, kind="stable")] = ordered
            if not np.array_equal(rows[:, 0], ops):
                raise ValueError("codec opcode plane disagrees")
            yield rows.tobytes()


def hash_file(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()
