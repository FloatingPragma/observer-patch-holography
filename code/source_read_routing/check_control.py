"""Small independent Python oracle, including every induced logical pair.

The production native verifier additionally enforces the frozen BFS route.
This oracle checks primitive semantics, the six-event protocol, the complete
read menu, and actual writer ancestry without using producer-assigned labels.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import struct

import numpy as np

from pack import segment_rows
from verify import check_specification, inspect_input, load, require

ABSENT = 2**64-1


def rows_for(packet, folder):
    for entry in packet["tape"]["segments"]:
        yield from (tuple(map(int, row)) for block in segment_rows(folder, entry["path"])
                    for row in np.frombuffer(block, dtype="<u8").reshape(-1, 8))


def check(packet, folder, rows=None):
    check_specification()
    q, n, rounds, centre, offsets, indices = inspect_input(folder/f"q{packet['inputs']['q']}.input", packet["inputs"])
    require(q == 3, "Python all-pairs oracle is the declared q3 control")
    variant = packet["variant"]
    c, ports = packet["inputs"]["carriers"], 12*packet["inputs"]["carriers"]
    raw = (folder/f"q{q}.input").read_bytes()
    host = np.frombuffer(raw, dtype="<u4", count=n, offset=64).tolist()
    z = np.frombuffer(raw, dtype="<i8", count=6*n, offset=len(raw)-48*n).reshape(n, 6).tolist()
    from verify import ROOT
    support = load(ROOT/packet["inputs"]["support_path"])
    seams = {tuple(sorted((12*a+p, 12*b+r))) for a,p,b,r in support["glued_pairs"]}
    prep = []
    site_at = dict(zip(host, range(n)))
    axes = [0,1,1,0,2,3,3,2,4,5,5,4]
    signs = [1,1,-1,-1,1,1,-1,-1,1,1,-1,-1]
    for cell in range(c):
        for port in range(12):
            value = 2*max(0, signs[port]*z[site_at[cell]][axes[port]]) if cell in site_at else 0
            prep.append((1, cell, value+(2*(port+1) if variant == "scratch" else 0)))
    prep.extend((0, cell, 0) for cell in range(c))
    accum = ports+c
    prep.extend((1, host[s], 0) for s in range(n))
    prep.extend((10, host[s], 2*z[s][a]) for s in range(n) for a in range(6))
    prep.extend((2, host[s], 2*(s+1+(variant == "source" and s == centre))) for s in range(n))
    state, ancestors, logical_masks, logical_values = [], [], [], []
    # State entries: scalar*2, writer, owner, immutable. Ancestry is calculated
    # from actual consumed writers, not from an asserted logical origin field.
    started = set()
    consumed = [set() for _ in range(n)]
    layer, committed, hop_phase, hops, scalar_reads, scalar_writes = 1, set(), 0, 0, 0, 0
    protocol = [3,4,5,6,4,4]
    active = None
    digest = hashlib.sha256()
    source = rows_for(packet, folder) if rows is None else rows
    maximum = 0
    for eid, row in enumerate(source):
        require(len(row) == 8, "event width")
        digest.update(struct.pack("<8Q", *row))
        kind,a,b,out,owner,wa,wb,encoded = row
        value = encoded if encoded < 2**63 else encoded-2**64
        mask = 0
        for key, writer in ((a,wa),(b,wb)):
            if key == ABSENT:
                require(writer == ABSENT, "extra writer")
                continue
            require(key < len(state) and writer == state[key][1] and writer < eid, "stale/absent writer")
            require(kind == 5 or state[key][2] == owner, "remote local read")
            mask |= ancestors[writer]
            scalar_reads += 1
        expected, immutable = None, kind in (0,2,6,9,10)
        if eid < len(prep):
            require((kind,owner,value) == prep[eid] and a == b == ABSENT and out == len(state), "prepared interface/value")
            expected = value
        else:
            require(layer <= rounds, "extra logical layer")
            if kind in (3,4,5,6):
                require(kind == protocol[hop_phase], "six-event transport protocol")
                if hop_phase == 0:
                    require(b == ABSENT and state[a][3] and out < ports, "export interface")
                    active = {"source": out, "archive": a}
                elif hop_phase == 1:
                    require(a == ports+owner and out < ports and b == ABSENT, "receiver reset interface")
                    active["receiver"] = out
                elif hop_phase == 2:
                    require(a == active["source"] and b == active["receiver"] and out == a and owner == ABSENT, "mean interface")
                    require(tuple(sorted((a,b))) in seams, "absent support seam")
                elif hop_phase == 3:
                    require(a == active["receiver"] and b == ABSENT and out == len(state), "capture interface")
                else:
                    endpoint = active["source"] if hop_phase == 4 else active["receiver"]
                    require(out == endpoint and owner == endpoint//12 and a == ports+owner and b == ABSENT, "cleanup interface")
                hop_phase = (hop_phase+1)%6
                if not hop_phase:
                    hops += 1
            else:
                require(hop_phase == 0, "interrupted hop")
            if kind == 3:
                expected = state[a][0]
            elif kind == 4:
                require(state[a][0] == 0 and state[a][3], "reset is not constant")
                expected = 0
            elif kind == 5:
                require(state[b][0] == 0 and (state[a][0]+state[b][0])%2 == 0, "mean domain")
                expected = (state[a][0]+state[b][0])//2
            elif kind == 6:
                expected = 2*state[a][0]
            elif kind == 7:
                s = out-accum
                require(0 <= s < n and s not in started and owner == host[s] and a == ports+owner and b == ABSENT, "accumulator start")
                started.add(s); expected = 2
            elif kind == 8:
                t = out-accum
                require(0 <= t < n and len(started) == n and a == out and owner == host[t] and state[b][3], "logical read interface")
                # Find the unique latest-layer logical writer contained in the
                # actual read's ancestry. Earlier ancestors are allowed.
                prior_mask = ancestors[state[b][1]] >> ((layer-1)*n)
                require(prior_mask and prior_mask & (prior_mask-1) == 0, "mixed/wrong-version transport")
                s = prior_mask.bit_length()-1
                require(s < n and s not in consumed[t], "duplicate read")
                consumed[t].add(s); expected = state[a][0]+state[b][0]
            elif kind == 9:
                t = a-accum
                require(0 <= t < n and t == len(committed) and owner == host[t] and b == ABSENT and out == len(state), "commit interface")
                require(consumed[t] == set(indices[offsets[t]:offsets[t+1]].tolist()), "missing/extra metric read")
                expected = state[a][0]+(2 if variant == "branch" and layer == 1 and t == centre else 0)
                committed.add(t)
                if len(committed) == n:
                    layer += 1; committed=set(); started=set(); consumed=[set() for _ in range(n)]
            else:
                raise ValueError("undeclared primitive")
        require(expected == value, "incorrect scalar value")
        if kind in (2,9):
            mask |= 1 << len(logical_masks)
            logical_masks.append(mask); logical_values.append(value//2)
        ancestors.append(mask)
        targets = (a,b) if kind == 5 else (out,)
        for target in targets:
            who = target//12 if kind == 5 else owner
            if target == len(state):
                state.append((value,eid,who,immutable))
            else:
                require(target < len(state) and not state[target][3], "immutable overwrite")
                require(state[target][2] == who, "remote write")
                state[target] = (value,eid,who,immutable)
            scalar_writes += 1
        maximum = max(maximum, abs(value))
    require(layer == rounds+1 and not hop_phase, "truncated execution")
    expected_masks = [1 << s for s in range(n)]
    expected_values = [s+1+(variant == "source" and s == centre) for s in range(n)]
    for j in range(1, rounds+1):
        for s in range(n):
            parents = indices[offsets[s]:offsets[s+1]].tolist()
            mask = 1 << (j*n+s)
            for p in parents:
                mask |= expected_masks[(j-1)*n+p]
            expected_masks.append(mask)
            expected_values.append(1+sum(expected_values[(j-1)*n+p] for p in parents)+(variant == "branch" and j == 1 and s == centre))
    require(logical_masks == expected_masks, "full induced logical order")
    require(logical_values == expected_values, "logical intervention values")
    costs = packet["costs"]
    for key, val in {"events": len(ancestors), "hops": hops, "registers": len(state),
                     "register_reads": scalar_reads, "register_writes": scalar_writes,
                     "max_abs_scaled_scalar": maximum}.items():
        require(costs[key] == val, "resource count: "+key)
    if rows is None:
        require(digest.hexdigest() == packet["tape"]["decoded_sha256"], "whole tape custody")
    return {"events": len(ancestors), "logical_pairs": len(logical_masks)**2, "logical_values": logical_values}
