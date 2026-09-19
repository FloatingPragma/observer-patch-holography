"""Independent support, exact metric, event, and intervention replay.

No producer or preparation module is imported. Native replay checks every
retained primitive and every writer; Python independently reconstructs the
complete metric relation and compares all completed logical values.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import cmp_to_key
import hashlib
import json
from math import isqrt
import os
from pathlib import Path
import struct
import subprocess
import tempfile

import numpy as np

from tape import decode, hash_file
from pack import segment_rows

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result
    def reject(value):
        raise ValueError("nonintegral JSON number: " + value)
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_float=reject, parse_constant=reject)


def check_specification(packet=None):
    packet = load(HERE/"specification.json") if packet is None else packet
    require(packet["schema"] == "oph.source_read_routing.specification.v1" and packet["issue"] == 777, "specification identity")
    require(packet["production_levels"] == [13,21] and packet["production_variants"] == ["baseline","source"], "complete production contract")
    require(packet["control_level"] == 3 and packet["control_variants"] == ["baseline","source","branch","scratch"], "control contract")
    expected = {"full_q13_q21_routed_execution":True, "axiomatic_read_law_derived":False,
                "universal_routing_impossibility":False, "physical_clock_identified":False,
                "quantum_instrument_implemented":False, "routing_events_counted_as_spacetime_volume":False}
    require(packet["scope"] == expected and all(type(x) is bool for x in packet["scope"].values()), "scope promotion")
    # Merely retaining the M1 label/owners is insufficient: the law and the
    # outstanding derivation must remain stated, not disappear at closeout.
    for key in ("population", "read_law", "local_feedback", "compiler", "derivation_obligation"):
        value = packet["M1"].get(key)
        require(isinstance(value, str) and bool(value.strip()), "missing M1 declaration: "+key)
    require(packet["M1"]["status"] == "retained supplied structural rule" and packet["M1"]["derivation_owners"] == [740,779]
            and packet["M1"]["premise_register"] == "PR-52", "M1 obligation transfer")
    require(packet.get("exit") == {
        "selected_route":"declared_M1_with_explicit_derivation_transfer",
        "acceptance_matrix":"code/source_read_routing/README.md#objective-deliverables-exit-and-audit-corrections"},
        "exit must retain M1 and its explicit derivation transfer")


def sign_sqrt5(a, b):
    # Scalar arbitrary-precision sign, independent of the producer's Q(phi) helper.
    if a*b >= 0:
        return (a+b > 0) - (a+b < 0)
    square = a*a - 5*b*b
    return ((a > 0) - (a < 0)) * ((square > 0) - (square < 0))


def refine(faces):
    next_vertex = max(max(face) for face in faces) + 1
    midpoint = {}
    out = []
    for a, b, c in faces:
        ids = []
        for u, v in ((a, b), (b, c), (c, a)):
            key = tuple(sorted((u, v)))
            if key not in midpoint:
                midpoint[key] = next_vertex
                next_vertex += 1
            ids.append(midpoint[key])
        ab, bc, ca = ids
        out.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]])
    return out


def check_support(packet, level):
    require(type(level) is int and level in (3, 4, 5), "support level")
    c = 20 * 4**level
    require(packet["level"] == level and packet["carriers"] == c, "support size")
    original = load(ROOT / "code/source_routing/support_w12_l3.json")
    expected = original["faces"]
    for _ in range(3, level):
        expected = refine(expected)
    require(packet["faces"] == expected, "not the committed spherical refinement")
    require(packet["intra_carrier_seams"] == original["intra_carrier_seams"], "local carrier changed")
    by_vertex = defaultdict(list)
    for cell, face in enumerate(expected):
        for vertex in face:
            by_vertex[vertex].append(cell)
    require(Counter(map(len, by_vertex.values())) == Counter({5: 12, 6: 10*4**level-10}), "spherical valences")
    wanted = {tuple(sorted((a, b))) for cells in by_vertex.values() for a in cells for b in cells if a != b}
    seen, ports = set(), set()
    for a, p, b, r in packet["glued_pairs"]:
        require(all(type(x) is int for x in (a, p, b, r)), "noninteger seam")
        require(0 <= a < b < c and 0 <= p < 12 and 0 <= r < 12, "seam bounds")
        require((a, b) not in seen and (a, p) not in ports and (b, r) not in ports, "duplicate seam/port")
        seen.add((a, b)); ports.update(((a, p), (b, r)))
    require(seen == wanted, "incomplete W12 shared-vertex support")
    require(len(seen) == 6*c-30, "glued seam census")


def inspect_input(path, info):
    data = Path(path).read_bytes()
    require(hashlib.sha256(data).hexdigest() == info["input_sha256"], "input digest")
    require(len(data) >= 64, "input header")
    q, n, c, rounds, centre, edges, reads, version = struct.unpack_from("<8Q", data)
    require(q in (3, 13, 21) and version == 1 and n == q**3, "input dimensions")
    level = {3: 3, 13: 4, 21: 5}[q]
    require(c == 20*4**level and rounds == isqrt(q)+(isqrt(q)**2 < q), "input support/rounds")
    offset = 64
    def take(dtype, count):
        nonlocal offset
        size = np.dtype(dtype).itemsize * count
        require(offset+size <= len(data), "truncated input array")
        result = np.frombuffer(data, dtype=dtype, count=count, offset=offset)
        offset += size
        return result
    host = take("<u4", n)
    pairs = take("<u4", edges*4).reshape(-1, 4)
    indptr = take("<u8", n+1)
    indices = take("<u4", reads)
    records = take("<i8", n*6).reshape(-1, 6)
    require(offset == len(data), "extra input bytes")
    require(indptr[0] == 0 and indptr[-1] == reads and np.all(indptr[1:] > indptr[:-1]) and np.all(indices < n), "metric CSR")
    support_relative = "code/source_routing/support_w12_l3.json" if level == 3 else f"code/source_read_routing/support_w12_l{level}.json"
    require(info["support_path"] == support_relative, "canonical support path")
    support_path = ROOT / support_relative
    require(hash_file(support_path) == info["support_sha256"], "support digest")
    support = load(support_path)
    check_support(support, level)
    require(np.array_equal(pairs, support["glued_pairs"]), "changed physical port assignment")
    coords = np.array([(a, b, d) for a in range(q) for b in range(q) for d in range(q)], dtype=np.int64)
    floor = np.array([(i+isqrt(5*i*i))//2 for i in range(q)], dtype=np.int64)
    m = -floor[coords]
    expected_records = np.column_stack((coords[:, 1]-m[:, 0], coords[:, 1]+m[:, 0],
                                        coords[:, 2]-m[:, 1], coords[:, 2]+m[:, 1],
                                        coords[:, 0]-m[:, 2], coords[:, 0]+m[:, 2]))
    require(np.array_equal(records, expected_records), "golden address records")
    order = sorted(range(q), key=cmp_to_key(lambda a,b: sign_sqrt5(
        int(-2*floor[a]+2*floor[b]+a-b), a-b)))
    rank = {b: a for a, b in enumerate(order)}
    keys = [sum(((rank[int(row[axis])] >> bit)&1) * 2**(3*bit+axis)
                for bit in range(q.bit_length()) for axis in range(3)) for row in coords]
    expected_host = np.empty(n, dtype=np.uint32)
    expected_host[np.argsort(keys, kind="stable")] = np.arange(n)
    require(np.array_equal(host, expected_host), "fixed injective Morton host assignment")
    distances = [(int(1+4*floor[b]**2+4*b*b+4*floor[b]), int(4*b*b-8*floor[b]*b-4*b)) for b in range(q)]
    best = min(range(q), key=cmp_to_key(lambda a,b: sign_sqrt5(
        2*(distances[a][0]-distances[b][0])+distances[a][1]-distances[b][1], distances[a][1]-distances[b][1])))
    require(centre == best*(q*q+q+1), "source intervention site")
    digest = hashlib.sha256(b"[")
    # All n^2 decisions, including equality/self, independently rebuilt in strips.
    for start in range(0, n, 32):
        da = m[start:start+32, None, :] - m[None, :, :]
        db = coords[start:start+32, None, :] - coords[None, :, :]
        constant = q*(2*(da*da+db*db).sum(axis=2)+(2*da*db+db*db).sum(axis=2))-2
        radical = q*(2*da*db+db*db).sum(axis=2)
        sign = np.where(constant*radical >= 0, np.sign(constant+radical),
                        np.sign(constant)*np.sign(constant*constant-5*radical*radical))
        for local, mask in enumerate(sign <= 0):
            s = start+local
            expected = np.flatnonzero(mask)
            require(np.array_equal(indices[int(indptr[s]):int(indptr[s+1])], expected), "metric read menu differs")
            if s:
                digest.update(b",")
            digest.update(("["+",".join(map(str, expected.tolist()))+"]").encode("ascii"))
    digest.update(b"]\n")
    require(digest.hexdigest() == info["neighbour_sha256"], "neighbour digest")
    for key, expected in {"q": q, "sites": n, "carriers": c, "rounds": rounds, "centre_site": centre, "logical_reads": reads*rounds}.items():
        require(type(info[key]) is int and info[key] == expected, "input census: "+key)
    if q in (13, 21):
        parent = json.loads((ROOT / "evidence/source_net_causal_poset/carrier_source_net_receipt.json").read_text(encoding="utf-8"))
        row = next(r for r in parent["levels"] if r["q"] == q)
        require(row["provenance"]["read_relation_sha256"] == digest.hexdigest() and row["operation_costs"]["total_reads"] == reads*rounds, "pinned family comparison")
    return q, n, rounds, centre, indptr.astype(np.int64), indices


def path_for_native(path):
    path = Path(path).resolve()
    return "/mnt/"+path.drive[0].lower()+path.as_posix()[2:] if os.name == "nt" else str(path)


def check_resources(packet):
    """Require the complete census before comparing it with native replay."""
    info = packet["inputs"]
    n, c, rounds = info["sites"], info["carriers"], info["rounds"]
    costs = packet["costs"]
    keys = {"events", "hops", "register_reads", "register_writes", "registers",
            "max_abs_scaled_scalar", "max_read_depth", "total_read_depth", "tree_hops_per_round"}
    require(set(costs) == keys and all(type(x) is int and x >= 0 for x in costs.values()), "incomplete semantic resource census")
    hops, reads = costs["hops"], info["logical_reads"]
    formulas = {"events":13*c+8*n+6*hops+reads+2*rounds*n,
                "register_reads":7*hops+2*reads+2*rounds*n,
                "register_writes":13*c+8*n+7*hops+reads+2*rounds*n,
                "registers":13*c+(rounds+8)*n+hops}
    require(all(costs[key] == value for key,value in formulas.items())
            and hops == rounds*costs["tree_hops_per_round"], "semantic resource formula")
    resources = packet["implementation_resources"]
    registers = costs["registers"]
    # Independently reconstruct the instrumented GNU/libstdc++ allocation
    # census; actual RSS is a retained host observation, not a deterministic
    # mathematical quantity to reproduce on another machine.
    capacity = 1 << (registers-1).bit_length()
    expected_resources = {"producer_cell_struct_bytes":24, "wide_arithmetic_bytes":16,
        "event_row_bytes":64, "codec_block_events":65536, "codec_raw_block_bytes":4194304,
        "register_heap_peak_requested_bytes":(capacity+capacity//2)*24,
        "auxiliary_heap_requested_bytes":72+76*n+252*c+
            ((c+63)//64)*8+16*(6*c-30)+4*(reads//rounds)}
    for key,value in expected_resources.items():
        require(type(resources.get(key)) is int and resources[key] == value, "implementation resource count: "+key)
    require(type(resources.get("native_peak_resident_bytes")) is int and resources["native_peak_resident_bytes"] > 0, "missing host memory observation")
    expected_work = {"bfs_vertex_dequeues":rounds*n*c,
        "bfs_directed_edge_examinations":rounds*n*2*(6*c-30),
        "bfs_nonroot_discoveries":rounds*n*(c-1),
        "parent_and_mark_initializations":2*rounds*n*c,
        "pruning_parent_traversals":packet["costs"]["hops"],
        "recipient_path_queries":info["logical_reads"],
        "tree_mark_examinations":rounds*n*(c-1)}
    require(packet.get("routing_control_work") == expected_work, "routing control work omitted or undercounted")


def verify(receipt_path, binary, work):
    check_specification()
    packet = load(receipt_path)
    require(packet["schema"] == "oph.source_read_routing.run.v1", "receipt schema")
    require(packet["producer_sha256"] == hash_file(HERE/"produce.cpp"), "producer source pin")
    variant = packet["variant"]
    require(variant in ("baseline", "source", "branch", "scratch"), "variant")
    folder = Path(receipt_path).parent
    info = packet["inputs"]
    input_path = folder / f"q{info['q']}.input"
    q, n, rounds, centre, indptr, indices = inspect_input(input_path, info)
    check_resources(packet)
    def blocks():
        if "segments" in packet["tape"]:
            segments = packet["tape"]["segments"]
            require(len(segments) == rounds+1, "complete layer segments")
            for phase, segment in enumerate(segments):
                name = f"q{q}_{variant}_phase{phase}.segment.json"
                require(segment["path"] == name and hash_file(folder/name) == segment["sha256"], "segment manifest custody")
                yield from segment_rows(folder, name)
        else:
            require(packet["tape"]["path"] == f"q{q}_{variant}.tape", "canonical tape path")
            tape_path = folder / packet["tape"]["path"]
            require(hash_file(tape_path) == packet["tape"]["sha256"] and tape_path.stat().st_size == packet["tape"]["bytes"], "tape custody")
            with tape_path.open("rb") as stream:
                yield from decode(stream)
    work.mkdir(parents=True, exist_ok=True)
    logical_path = work / f"q{q}_{variant}.logical"
    diagnostics = work / f"q{q}_{variant}.verify.stderr"
    prefix = ["wsl", "-d", "Ubuntu", "--"] if os.name == "nt" else []
    cmd = prefix + [path_for_native(binary), path_for_native(input_path), variant, path_for_native(logical_path)]
    digest, size = hashlib.sha256(), 0
    with diagnostics.open("wb") as error:
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=error)
        try:
            for data in blocks():
                digest.update(data); size += len(data)
                process.stdin.write(data)
            process.stdin.close()
            result = process.stdout.read()
        except (BrokenPipeError, OSError):
            process.wait()
            raise ValueError(diagnostics.read_text(encoding="utf-8")) from None
        except BaseException:
            process.kill(); process.wait(); raise
        require(process.wait() == 0, diagnostics.read_text(encoding="utf-8"))
    require(size == packet["tape"]["decoded_bytes"] and size//64 == packet["tape"]["events"] and digest.hexdigest() == packet["tape"]["decoded_sha256"], "decoded tape custody")
    actual = json.loads(result)
    for key, value in packet["costs"].items():
        require(actual[key] == value, "resource undercount: "+key)
    require(actual["logical_reads"] == info["logical_reads"], "incomplete logical family")
    logical = np.fromfile(logical_path, dtype="<i8")
    require(logical.size == n*(rounds+1), "logical output size")
    logical = logical.reshape(rounds+1, n)
    expected = np.arange(1, n+1, dtype=np.int64)
    if variant == "source":
        expected[centre] += 1
    require(np.array_equal(logical[0], expected), "logical preparation")
    for layer in range(1, rounds+1):
        # For these declared bounded levels the exact positive recurrence fits int64.
        degree = int(np.max(np.diff(indptr)))
        require(int(expected.max())*degree+2 < 2**63, "logical reference overflow")
        expected = 1 + np.add.reduceat(expected[indices], indptr[:-1])
        if variant == "branch" and layer == 1:
            expected[centre] += 1
        require(np.array_equal(logical[layer], expected), "intervention-preserving logical refinement")
    actual.update({"logical_value_sha256": hash_file(logical_path), "metric_pairs_checked": n*n,
                   "decoded_sha256": digest.hexdigest()})
    print(json.dumps(actual, sort_keys=True), flush=True)
    return actual


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    if args.work:
        verify(args.receipt, args.binary, args.work)
    else:
        with tempfile.TemporaryDirectory(prefix="oph-routing-verify-") as directory:
            verify(args.receipt, args.binary, Path(directory))
