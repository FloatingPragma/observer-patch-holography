"""Compile the eight phase contexts and their unsharp state-update instruments."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import fixed_pointer, sequential
from ..full_address import measure as address

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
TABLE = ROOT/"code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json"


def context_table():
    return json.loads(TABLE.read_text())["contexts"]


def projector(context):
    row = next(row for row in context_table() if row["context"] == context)
    return np.array([[float(Fraction(v[0]))+float(Fraction(v[1]))*np.sqrt(3)
                      +1j*(float(Fraction(v[2]))+float(Fraction(v[3]))*np.sqrt(3))
                      for v in line] for line in row["effect"]])


def pulse(kind,a,b,ratio):
    return kind,a,b,Fraction(ratio)


def path(a,b):
    board = federation(0)
    via = {a:None}
    queue = [a]
    for current in queue:
        if b in via:
            break
        for target in board.neighbours[current]:
            if target<12 and target not in via:
                via[target] = current
                queue.append(target)
    sequential.require(b in via,"data route absent")
    result = [b]
    while result[-1] != a:
        result.append(via[result[-1]])
    return list(reversed(result))


def routed_contact(a,b,ratio):
    route = path(a,b)
    transport = [("contact",x,y,Fraction(1))
                 for x,y in reversed(list(zip(route[1:-1],route[2:])))]
    return transport+[pulse("contact",a,route[1],ratio)]+list(reversed(transport))


def transposition(a,b):
    route = path(a,b)
    pairs = list(zip(route,route[1:]))
    return [("contact",x,y,Fraction(1)) for x,y in pairs+list(reversed(pairs[:-1]))]


def context_word(context):
    sequential.require(context in {r["context"] for r in context_table()},"unknown phase context")
    if context in ("web_diagonal","web_conjugated_0","web_conjugated_1"):
        return []
    if context == "phase":
        return [pulse("contact",0,1,Fraction(3,2))]
    dwell = Fraction(4,3) if context in ("web_conjugated_2","web_conjugated_3") else Fraction(2,3)
    return [pulse("reference",0,0,Fraction(1,16)),pulse("contact",0,1,dwell),
            pulse("reference",0,0,Fraction(3,16))]


def normalize_mode(mode):
    return "coarse" if mode == "weak" else mode


def dilation_word(mode):
    mode = normalize_mode(mode)
    sequential.require(mode in ("sharp","coarse","fine","coin"),"unknown instrument mode")
    waits = {"sharp":(0,1),"coarse":(Fraction(1,3),Fraction(2,3)),
             "fine":(Fraction(1,3),Fraction(2,3)),"coin":(Fraction(1,3),Fraction(1,3))}[mode]
    return [step for j,dwell in enumerate(waits) if dwell
            for step in routed_contact(j,j+2,dwell)]


def keep_mask(mode):
    # Keep high bits as explicit out-of-dilation rejection records. Only the
    # logical bit zero is uncomputed in the coarse instrument.
    return 15 if normalize_mode(mode)=="fine" else 14


def record_word(mode):
    calls = address.schedule()
    unread = [call for call in calls if not (keep_mask(mode)>>call["bit"])&1]
    return address.pulse_word(calls)+fixed_pointer.inverse(address.pulse_word(unread))


def post_word(context,mode,raw):
    mode = normalize_mode(mode)
    accepted = (0,1,2,3) if mode=="fine" else (0,2)
    if raw not in accepted:
        return []  # Complete rejection branches retain their actual source.
    y = (raw>>1)&1
    decode = transposition(0,2)+transposition(1,3) if y else []
    correction = fixed_pointer.diagonal_word(0,7,Fraction(5,51)) if mode in ("coarse","fine") else []
    return decode+correction+fixed_pointer.inverse(context_word(context))


@lru_cache(maxsize=None)
def source_unitary(word,sign):
    """Actual one-defect columns for the same native source word on both signs."""
    sequential.require(sign in (-1,1),"signed source required")
    board = federation(0)
    result = np.eye(12,dtype=complex)
    for kind,a,b,ratio in word:
        sequential.require(0<=a<12 and 0<=b<12,"source pulse outside data cell")
        if kind=="contact":
            sequential.require(b in board.neighbours[a],"source contact outside incidence")
            h = np.zeros((12,12))
            h[np.ix_([a,b],[a,b])] = [[8.5,-0.5],[-0.5,8.5]]
        else:
            sequential.require(kind=="reference" and board.degree[a]==5,"native source reference absent")
            h = 8*np.eye(12)
            h[a,a] = 8*(1+sign)**2
        result = expm(-1j*float(ratio)*np.pi*h)@result
    return result


def raw_branch_maps(context,mode="sharp",sign=-1):
    """Seventeen physical branches H_2 -> H_12, with complement retained."""
    mode = normalize_mode(mode)
    prefix = context_word(context)+dilation_word(mode)
    columns = source_unitary(tuple(prefix),sign)[:,:2]
    kept = sum(1 for call in address.schedule() if (keep_mask(mode)>>call["bit"])&1)
    scalar = address.COMMON_PHASE**kept
    branches = {}
    for raw in range(17):
        selected = columns.copy()
        selected[[a&keep_mask(mode)!=raw for a in range(12)]] = 0
        branches[raw] = scalar*source_unitary(tuple(post_word(context,mode,raw)),sign)@selected
    return branches


def branch_matrices(context,mode="sharp",sign=-1):
    """Native yes/no H_2 -> H_2 matrices for sharp, weak/coarse or coin.

    sign selects an audit block of one common physical word; it is not a
    controller input. Fine mode has two Kraus branches per grouped label.
    """
    mode = normalize_mode(mode)
    sequential.require(mode in ("sharp","coarse","coin"),"use raw_branch_maps for the fine instrument")
    raw = raw_branch_maps(context,mode,sign)
    sequential.require(max(np.linalg.norm(k[2:]) for k in raw.values())<3e-12,
                       "conditional decode left the two logical rails")
    return {0:raw[0][:2],1:raw[2][:2]}


def fine_matrices(context,sign=-1):
    return {raw:k[:2] for raw,k in raw_branch_maps(context,"fine",sign).items() if raw<4}


def effect(context,mode="sharp"):
    mode = normalize_mode(mode)
    p = projector(context)
    return p if mode=="sharp" else (3*np.eye(2)/4 if mode=="coin" else np.eye(2)/4+p/2)
