"""Reexecute, prove byte identity, and attach actual host allocation accounting.

GNU time's peak resident set is a host-software measurement. It is not a
physical observer capacity or an empirical physics result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from run import native_path
from tape import hash_file

HERE = Path(__file__).resolve().parent


def account(receipt_path, binary, work):
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    q, variant = receipt["inputs"]["q"], receipt["variant"]
    work.mkdir(parents=True, exist_ok=True)
    error_path, memory_path = work/f"q{q}_{variant}.account.stderr", work/f"q{q}_{variant}.rss"
    input_path = receipt_path.parent/f"q{q}.input"
    # The native event stream is still produced in full; only its second
    # physical copy is replaced by a byte-for-byte SHA-256 comparison.
    import os
    prefix = ["wsl", "-d", "Ubuntu", "--"] if os.name == "nt" else []
    args = prefix+["/usr/bin/time", "-f", "%M", "-o", native_path(memory_path),
                   native_path(binary), native_path(input_path), variant]
    digest, count = hashlib.sha256(), 0
    with error_path.open("wb") as errors:
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=errors)
        while data := process.stdout.read(4*1024*1024):
            digest.update(data); count += len(data)
        if process.wait():
            raise ValueError(error_path.read_text(encoding="utf-8"))
    if digest.hexdigest() != receipt["tape"]["decoded_sha256"] or count != receipt["tape"]["decoded_bytes"]:
        raise ValueError("accounting rerun changed the retained execution")
    reported = json.loads(error_path.read_text(encoding="utf-8").splitlines()[-1])
    for key, value in receipt["costs"].items():
        if reported[key] != value:
            raise ValueError("accounting rerun changed semantic cost: "+key)
    resources = {key: value for key, value in reported.items() if key not in receipt["costs"]}
    resources.update({"native_peak_resident_bytes": int(memory_path.read_text().strip())*1024,
                      "native_peak_resident_scope": "GNU time %M on the executing Linux host; includes allocator, stack and runtime; not a physical patch capacity.",
                      "register_heap_scope": "Requested vector storage including simultaneous old/new allocation during growth; allocator bookkeeping excluded here and included in host peak RSS.",
                      "codec_block_events": 65536, "codec_raw_block_bytes": 4194304})
    receipt["implementation_resources"] = resources
    n, c, k = receipt["inputs"]["sites"], receipt["inputs"]["carriers"], receipt["inputs"]["rounds"]
    receipt["routing_control_work"] = {
        "bfs_vertex_dequeues": k*n*c,
        "bfs_directed_edge_examinations": k*n*2*(6*c-30),
        "bfs_nonroot_discoveries": k*n*(c-1),
        "parent_and_mark_initializations": 2*k*n*c,
        "pruning_parent_traversals": receipt["costs"]["hops"],
        "recipient_path_queries": receipt["inputs"]["logical_reads"],
        "tree_mark_examinations": k*n*(c-1)}
    receipt["producer_sha256"] = hash_file(HERE/"produce.cpp")
    receipt_path.write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"q":q,"variant":variant,"byte_identical":True,**resources},sort_keys=True),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    account(args.receipt,args.binary,args.work)
