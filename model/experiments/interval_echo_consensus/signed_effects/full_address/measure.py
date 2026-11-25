"""Encode all twelve defect addresses in four fixed material pointer bits."""

from __future__ import annotations

import argparse
from collections import deque
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import fixed_pointer, sequential

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
RECEIPT = HERE / "results.json"
SERVICE_CELLS = (1, 2, 4, 5)
REFERENCE = 239
SITES = tuple({"source": 12*c+1, "park": 12*c+4,
               "cross": 12*c+5, "host": 12*c+7} for c in SERVICE_CELLS)
RESERVED = frozenset({REFERENCE, *(site[role] for site in SITES
                                 for role in ("park", "cross", "host"))})
COMMON_PHASE = np.exp(-133j*np.pi/260)


def route(board, start, end):
    sequential.require(start not in RESERVED and end not in RESERVED,
                       "route endpoint is reserved")
    parents = {start: None}
    queue = deque([start])
    while end not in parents:
        sequential.require(bool(queue), "no isolated native route")
        previous = queue.popleft()
        for target in board.neighbours[previous]:
            if target not in RESERVED and target not in parents:
                parents[target] = previous
                queue.append(target)
    path = [end]
    while path[-1] != start:
        path.append(parents[path[-1]])
    return list(reversed(path))


def local_word(site):
    aliases = {0: site["source"], 6: site["park"],
               10: site["cross"], 7: site["host"]}
    return [(kind, aliases[a], aliases[b], duration)
            for kind,a,b,duration in fixed_pointer.copy_word(0)]


def schedule(control="native"):
    board = federation(0)
    calls = []
    for bit, site in enumerate(SITES):
        sequential.require(board.degree[site["source"]] == 5
                           and board.degree[site["cross"]] == 5,
                           "native degree reference unavailable")
        for address in range(12):
            if not (address >> bit) & 1 or (control == "omit_bit3" and bit == 3):
                continue
            actual = 9 if control == "wrong_route" and bit == 3 and address == 8 else address
            path = route(board, actual, site["source"])
            for kind,a,b,duration in local_word(site):
                if kind == "contact":
                    sequential.require(b in board.neighbours[a], "non-native service contact")
                else:
                    sequential.require(board.degree[a] == 5, "degree-six reference")
            calls.append({"bit": bit, "declared_address": address,
                          "queried_address": actual, "site": bit, "route": path})
    return calls


def pulse_word(calls):
    result = []
    for call in calls:
        edges = list(zip(call["route"],call["route"][1:]))
        result += [("contact",a,b,fixed_pointer.Fraction(1)) for a,b in edges]
        result += local_word(SITES[call["site"]])
        result += [("contact",a,b,fixed_pointer.Fraction(1)) for a,b in reversed(edges)]
    return result


@lru_cache(maxsize=1)
def local_check():
    """Full 3^4 operator check; no projection during the native local word."""
    law = candidate.native_law
    h0,h1 = law.local_term_parts()
    full = h0 + candidate.KAPPA*h1
    h = full[:9,:9]
    sequential.require(np.linalg.norm(full[9:,:9],2) < 4e-13,
                       "native blank-memory sector not reducing")
    swap = np.eye(9)[:,[law.SRC.index((b,a)) for a,b in law.SRC]]
    swap_error = float(np.linalg.norm(expm(-1j*np.pi*full)[:,:9]
                                      - np.vstack((swap,np.zeros((90,9)))),2))
    columns = np.eye(81,dtype=complex).reshape(3,3,3,3,81)
    roles = {0:0,6:1,10:2,7:3}
    for kind,a,b,duration in fixed_pointer.copy_word(0):
        axes = [roles[a]] if kind == "reference" else [roles[a],roles[b]]
        local = np.diag([0,8,32]) if kind == "reference" else h
        gate = expm(-1j*float(duration)*np.pi*local)
        moved = np.moveaxis(columns,axes,range(len(axes)))
        updated = gate @ moved.reshape(len(gate),-1)
        columns = np.moveaxis(updated.reshape(moved.shape),range(len(axes)),axes)
    unitary = columns.reshape(81,81)
    encode = np.zeros((81,6),complex)
    expected = np.zeros_like(encode)
    for load in range(3):
        for bit in range(2):
            pointer = (2,0) if bit == 0 else (0,2)
            output = bit ^ int(load != 1)
            changed = (2,0) if output == 0 else (0,2)
            encode[np.ravel_multi_index((load,*pointer,1),(3,)*4),2*load+bit] = 1
            expected[np.ravel_multi_index((load,*changed,1),(3,)*4),2*load+bit] = COMMON_PHASE
    error = float(np.linalg.norm(unitary@encode-expected,2))
    sequential.require(error < 2e-12 and swap_error < 2e-12,
                       "local native material-copy or ternary swap identity failed")
    return {"full_local_dimension":81,"source_pointer_code_columns":6,
            "local_copy_operator_error":error,"native_pi_swap_operator_error":swap_error}


def initial_word(sign,address,pointer_word=0):
    loads = np.ones(240,dtype=int)
    loads[address] = 1+sign
    loads[REFERENCE] = 1-sign
    for bit,site in enumerate(SITES):
        loads[site["park"]],loads[site["cross"]] = ((0,2) if (pointer_word>>bit)&1 else (2,0))
    return loads


def execute_basis(calls,sign,address,pointer_word=0):
    """Exact operator-factorization execution on a physical board word.

    Every routed permutation is executed. The local code identity is checked
    on all six columns in local_check(), including one common phase and all
    81 output rows; its analytic factorization permits this compression.
    """
    loads = initial_word(sign,address,pointer_word)
    for call in calls:
        path = call["route"]
        for a,b in zip(path,path[1:]):
            loads[a],loads[b] = loads[b],loads[a]
        site = SITES[call["site"]]
        sequential.require(loads[site["host"]] == 1
                           and (loads[site["park"]],loads[site["cross"]]) in ((2,0),(0,2)),
                           "route disturbed an ancillary preparation")
        if loads[site["source"]] != 1:
            a,b = site["park"],site["cross"]
            loads[a],loads[b] = loads[b],loads[a]
        for a,b in reversed(list(zip(path,path[1:]))):
            loads[a],loads[b] = loads[b],loads[a]
    return loads,COMMON_PHASE**len(calls)


def read_word(loads):
    bits = []
    for site in SITES:
        pair = tuple(loads[[site["park"],site["cross"]]])
        if loads[site["host"]] != 1 or pair not in ((2,0),(0,2)):
            return 16
        bits.append(int(pair == (0,2)))
    return sum(bit<<index for index,bit in enumerate(bits))


def encoder(calls):
    """Exact reduced isometry on 24 input states and all 17 output labels."""
    local_check()
    v = np.zeros((24,17,24),complex)
    all_restored = True
    for sign_index,sign in enumerate((-1,1)):
        for address in range(12):
            for prior in range(16):
                loads,phase = execute_basis(calls,sign,address,prior)
                record = read_word(loads)
                expected = initial_word(sign,address,record)
                all_restored &= bool(np.array_equal(loads,expected))
                if prior == 0:
                    col = 12*sign_index+address
                    v[col,record,col] = phase
    sequential.require(all_restored,"data or ancillary route restoration failed")
    return v,all_restored


def sources():
    names = ["model/candidate_laws/contrast_memory_spherical_v1/law.py",
             "model/candidate_laws/interval_echo/candidate.py",
             "model/executable_model/carrier.py",
             "model/executable_model/carriers/FEDERATION_LEVEL0.npz",
             "model/experiments/interval_echo_consensus/signed_effects/fixed_pointer.py",
             "model/experiments/interval_echo_consensus/signed_effects/full_address/measure.py"]
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def report():
    check = local_check()
    board = federation(0)
    calls = schedule()
    word = pulse_word(calls)
    v,restored = encoder(calls)
    ideal = np.zeros_like(v)
    for index in range(24):
        ideal[index,index%12,index] = COMMON_PHASE**20
    error = float(np.linalg.norm((v-ideal).reshape(24*17,24),2))
    sequential.require(error < 2e-12,"full address operator map failed")
    controls = {}
    for name in ("omit_bit3","wrong_route"):
        wrong,_ = encoder(schedule(name))
        controls[name] = {"source_address":8,
                          "vacancy_records":np.sum(abs(wrong[:,:,8])**2,axis=0).tolist(),
                          "excess_records":np.sum(abs(wrong[:,:,20])**2,axis=0).tolist()}
        sequential.require(controls[name]["vacancy_records"][0] > 1-3e-12
                           and controls[name]["excess_records"][0] > 1-3e-12,
                           "address control missed its wrong label")
    # A native pi/2 contact prepares a two-address complex superposition.
    # Coherent inverse playback restores it; reading the pointer first pinches
    # those two branches and halves the inverse-analyzer return probability.
    source_contact = np.array([[8.5,-0.5],[-0.5,8.5]])
    prep = expm(-0.5j*np.pi*source_contact)[:,0]
    undo = expm(-1.5j*np.pi*source_contact)
    coherent_return = float(abs((undo@prep)[0])**2)
    dephased_return = float(sum(abs((undo[:,index]*prep[index])[0])**2 for index in range(2)))
    shots = 4096
    banks = []
    for sign_index,sign in enumerate((-1,1)):
        for address in range(12):
            counts = [0]*17
            counts[address] = shots
            banks.append({"sign":sign,"prepared_address":address,
                          "all_17_record_counts":counts,"shots":shots})
    contact_count = sum(kind=="contact" for kind,_,_,_ in word)
    reference_count = len(word)-contact_count
    return {
        "scope":"one unknown-sign defect at any of twelve data-cell addresses",
        "data_ports":list(range(12)),"balancing_reference_port":REFERENCE,
        "carrier_ports":board.n,"carrier_seams_sha256":hashlib.sha256(json.dumps(board.seams).encode()).hexdigest(),
        "service_sites":list(SITES),"reserved_ports":sorted(RESERVED),
        "all_service_reference_degrees":[int(board.degree[site[role]]) for site in SITES for role in ("source","cross")],
        "ordered_calls":calls,
        "native_pulse_word":[[kind,a,b,str(t)] for kind,a,b,t in word],
        "native_waits":len(word),"contact_waits":contact_count,"reference_waits":reference_count,
        "exact_reduction":"ternary pi SWAP, full81-state local code identity, routed conjugation, linearity",
        "literal_execution":"local matrix checks and all384 physical-board code words; no3^240 propagation or literal event count",
        "all_data_and_ancillary_ports_restored":restored,"local_checks":check,
        "reduced_operator_error":error,
        "numerical_telescope_diagnostic":20*check["local_copy_operator_error"]
            +(contact_count-20*10)*check["native_pi_swap_operator_error"],
        "source_sha256":sources(),
        "record_alphabet":"binary words0..15 and invalid16; invalid aggregates any malformed pointer or host",
        "unused_words":[12,13,14,15],"coherent_native_analyzer_return":coherent_return,
        "pointer_dephased_native_analyzer_return":dephased_return,
        "controls":controls,"banks":banks,
        "pulse_error_contract":f"complete cq half-trace error <= min(1,epsilon_prep+{33*contact_count+32*reference_count}*delta) for absolute duration error delta per wait",
        "final_read":"supplied complete sharp occupations of four pointer dimers and hosts; no intermediate read",
        "custody":"producer-side software binding only",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        with RECEIPT.open("x") as output:
            output.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        sequential.compare(json.loads(RECEIPT.read_text()),result)
        print("full-cell address register: all twelve labels and native routes replay")


if __name__ == "__main__":
    main()
