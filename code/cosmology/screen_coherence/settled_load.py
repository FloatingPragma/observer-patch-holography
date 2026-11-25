"""Measure load and contact covariance in the archived strict-repair observer records."""

import argparse
import hashlib
import json
import struct
import sys
from pathlib import Path

import numpy as np
from numba import njit

try:
    from isal.isal_zlib import decompress
except ImportError:
    from zlib import decompress

ARCHIVE = "https://dovybzrnzxfo1.cloudfront.net"
MANIFEST_SHA = "3c0b240ac4464afd9b82251b46a86c6622a87a8a281d0ad159b641413360f824"
FIELDS = ("initial_load", "settled_load", "change")
REGION_SHIFTS = (0, 6, 12, 20)


@njit
def add_pair(moments, a, b):
    for k in range(3):
        x, y = np.int64(a[k]), np.int64(b[k])
        moments[k, 0] += 1
        moments[k, 1] += x + y
        moments[k, 2] += x*x + y*y
        moments[k, 3] += x*y


@njit
def measure_shard(raw, first, neighbors):
    n = len(raw)
    values = np.empty((n, 3), np.int8)
    pairs = np.zeros((2, 3, 4), np.int64)
    scales = np.zeros((4, 7), np.int64)
    for i in range(n):
        mask = int(raw[i, 0]) + 256 * int(raw[i, 1])
        initial, change = 0, 0
        for p in range(12):
            bit = (mask >> p) & 1
            initial += 2 * bit
            if raw[i, 2+2*p] or raw[i, 3+2*p]:
                change += 1 - 2 * bit
        values[i] = (initial-12, initial+change-12, change)
    for i in range(n):
        owner = first + i
        control = i ^ (n // 2)
        touches_control = False
        for j in range(3):
            peer = np.int64(neighbors[i, j])
            if first <= peer < owner:
                add_pair(pairs[0], values[i], values[peer-first])
            if peer == first + control:
                touches_control = True
        if control < i and not touches_control:
            add_pair(pairs[1], values[i], values[control])
    # A region is the descendants of one native refinement face. Count its
    # actual boundary links and signed writes from the member reports.
    for k, shift in enumerate(REGION_SHIFTS):
        size = 1 << shift
        for start in range(0, n, size):
            a, d, boundary, flux = 0, 0, 0, 0
            for i in range(start, start+size):
                a += values[i, 0]
                d += values[i, 2]
                for j in range(3):
                    if (neighbors[i, j] >> shift) == ((first+i) >> shift):
                        continue
                    boundary += 1
                    p = raw[i, 44+j] & 15
                    mate = (raw[i, 26+p//2] >> (4*(p%2))) & 15
                    if mate == 12:
                        flux += 1 - 2*((raw[i, p//8] >> (p%8)) & 1)
            if d != flux or abs(d) > boundary:
                raise ValueError("Regional change differs from recorded boundary flux")
            scales[k, 0] += 1
            scales[k, 1] += a
            scales[k, 2] += d
            scales[k, 3] += a*a
            scales[k, 4] += (a+d)*(a+d)
            scales[k, 5] += d*d
            scales[k, 6] += boundary
    return values, pairs, scales


def records(path, entry, rows, frame_size):
    data = path.read_bytes()
    if len(data) != entry["bytes"] or hashlib.sha256(data).hexdigest() != entry["sha256"]:
        raise ValueError(f"Archive checksum failed for {path}")
    frames = rows // frame_size
    offsets = struct.unpack_from(f"<{frames+1}I", data)
    if offsets[0] != 4*(frames+1) or offsets[-1] != len(data):
        raise ValueError("Invalid frame table")
    raw = np.empty((rows, 48), np.uint8)
    for i, (a, b) in enumerate(zip(offsets, offsets[1:])):
        frame = np.frombuffer(decompress(data[a:b]), dtype=np.uint8)
        raw[i*frame_size:(i+1)*frame_size] = frame.reshape(frame_size, 48)
    return raw


def describe_pairs(moments):
    result = {}
    for field, (number, total, squares, products) in zip(FIELDS, moments.tolist()):
        mean = total / (2 * number)
        variance = squares / (2 * number) - mean**2
        covariance = products / number - mean**2
        result[field] = dict(pairs=number, covariance=covariance,
                             correlation=covariance/variance)
    return result


def analyze(folder):
    manifest_data = (folder / "manifest.json").read_bytes()
    if hashlib.sha256(manifest_data).hexdigest() != MANIFEST_SHA:
        raise ValueError("This analysis is registered to the published level-13 manifest")
    m = json.loads(manifest_data)
    n = m["observers_per_shard"]
    pairs = np.zeros((2, 3, 4), np.int64)
    scales = np.zeros((4, 7), np.int64)
    across = {}
    for index, entry in enumerate(m["shards"]):
        raw = records(folder / entry["file"], entry, n, m["observers_per_frame"])
        neighbors = raw[:, 32:44].view("<u4")
        first = index*n
        values, local_pairs, local_scales = measure_shard(raw, first, neighbors)
        pairs += local_pairs
        scales += local_scales
        # Retain only boundary endpoints until the reciprocal shard arrives.
        ii, jj = np.nonzero((neighbors < first) | (neighbors >= first+n))
        left, right = [], []
        for i, j in zip(ii.tolist(), jj.tolist()):
            owner, peer = first+i, int(neighbors[i,j])
            if owner < peer:
                across[(owner,peer)] = values[i].copy()
            else:
                left.append(across.pop((peer,owner)))
                right.append(values[i])
        if left:
            accumulate_boundary(pairs[0], np.array(left), np.array(right))
        if (index+1) % 64 == 0:
            print(f"Read {index+1} of {len(m['shards'])} shards", file=sys.stderr, flush=True)
    if across or pairs[0,0,0] != m["external_seams"]:
        raise ValueError("Incomplete census of contact pairs")
    regions = []
    for shift, row in zip(REGION_SHIFTS, scales.tolist()):
        count, initial, change, aa, ff, dd, boundary = row
        size = 1 << shift
        if count*size != m["observers"] or change or initial+12*m["observers"] != m["initial_and_final_load"]:
            raise ValueError("Regional records fail the global load check")
        centered = initial / count
        regions.append(dict(observers_per_region=size, regions=count,
            mean_boundary_links=boundary/count, mean_load_per_observer=12+centered/size,
            variance_initial_total_per_member=(aa/count-centered**2)/size,
            variance_settled_total_per_member=(ff/count-centered**2)/size,
            variance_change_total_per_member=dd/count/size))
    return dict(archive=ARCHIVE+"/manifest.json", manifest_sha256=MANIFEST_SHA,
                observers=m["observers"], seed=m["seed"], regions=regions,
                contact_pairs=describe_pairs(pairs[0]), noncontact_control=describe_pairs(pairs[1]),
                control="Within each shard pair opposite ID halves at the same offset, excluding direct contacts.",
                region_check="Every region's change equals its signed boundary transfers.")


@njit
def accumulate_boundary(moments, left, right):
    for i in range(len(left)):
        add_pair(moments, left[i], right[i])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.archive), indent=2, allow_nan=False))
