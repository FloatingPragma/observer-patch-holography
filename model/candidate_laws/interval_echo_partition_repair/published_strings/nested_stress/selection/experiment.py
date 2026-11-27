"""Archive every bath record for two correlated native memories and a reference."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

from ...memories import Pair, THETA
from ...native_pair import NativePair, endpoint, entropy
from .service import SWEEP, native_gate, output_without_conditioning, ready_records, replacement_horizon, stock_read

HERE = Path(__file__).resolve().parent


def trace_distance(a, b):
    delta = a-b
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta+delta.conj().T)/2)))/2)


def swap_memory_order(a):
    return a.reshape(4,4,2,4,4,2).transpose(1,0,2,4,3,5).reshape(32,32)


def reference(a):
    return np.trace(a.reshape(16,2,16,2), axis1=0, axis2=2)


def correlated_source():
    # Distinct nonseparable correlations to a two-word exterior reference.
    v = np.arange(1,33) + 1j*((17*np.arange(32)+3) % 29)
    v = v/np.linalg.norm(v)
    return np.outer(v, v.conj())


def main_banks():
    codes = (NativePair(0, "right"), NativePair(1, "left"))
    reads = [stock_read("held", shape, 32) for shape in ("rest", "xy")]
    works = [row["work"] for row in reads]
    targets = [code.state(THETA+Pair.from_weighted_work(w).angle_increment)
               for code,w in zip(codes,works,strict=True)]
    target = np.kron(*targets)
    data = correlated_source()
    ab = np.trace(data.reshape(4,4,2,4,4,2), axis1=2, axis2=5)
    partial_transpose = ab.transpose(0,3,2,1).reshape(16,16)
    entanglement_witness = float(np.min(np.linalg.eigvalsh(partial_transpose)))
    assert entanglement_witness < -1e-3
    records = list(product((0,1), repeat=3))
    probabilities = []
    total = np.zeros_like(data)
    worst = 0.
    for first_block in ready_records(data, works[0], codes[0]).values():
        first = swap_memory_order(first_block)
        row = []
        for second_block in ready_records(first, works[1], codes[1]).values():
            result = swap_memory_order(second_block)
            total += result
            proba = float(np.trace(result).real)
            row.append(round(proba, 12))
            if proba > 1e-14:
                distance = trace_distance(result, np.kron(target, reference(result)))/proba
                worst = max(worst, distance)
            else:
                assert np.linalg.norm(result) < 1e-13
        probabilities.append(row)
    assert worst < 4e-13
    assert trace_distance(total, np.kron(target, reference(data))) < 4e-13
    assert abs(np.trace(total)-1) < 4e-13
    gains = []
    for code, w, ready in zip(codes, works, targets, strict=True):
        gain = entropy(endpoint(ready, "left"))-entropy(endpoint(code.state(THETA), "left"))
        assert abs(gain-Pair.from_weighted_work(w).cut_entropy_change) < 3e-13
        gains.append(round(gain, 12))
    return dict(survey=[{k:round(v,12) if isinstance(v,float) else v for k,v in r.items()} for r in reads],
                visible_records=["".join(map(str,z)) for z in records], joint_record_probabilities=probabilities,
                input_memories_partial_transpose_min_eigenvalue=round(entanglement_witness,12),
                nonzero_record_count=sum(p>0 for row in probabilities for p in row),
                worst_conditional_trace_distance_less_than=4e-13,
                summed_reference_marginal_error_less_than=4e-13, selected_cut_entropy_gains=gains)


def forgotten_address_control():
    pair = NativePair(0, "right")
    gate = native_gate(pair, 0)
    states = [np.outer(gate[:,i], gate[:,i].conj()) for i in (0,3)]
    rows = []
    for sweep in (SWEEP, (0,1,0)):
        outputs = [output_without_conditioning(rho, sweep, 0, pair) for rho in states]
        probabilities = [np.diag(gate.conj().T @ rho @ gate).real for rho in outputs]
        rows.append(dict(addresses=list(sweep),
                         high_low_word_probabilities=[[max(0,round(float(x),12)) for x in row] for row in probabilities],
                         distance_between_preparations=round(trace_distance(*outputs),12)))
    assert rows[0]["distance_between_preparations"] == 0
    assert rows[1]["distance_between_preparations"] == 0.75
    return rows


def deposit():
    horizons = []
    for levels in (1, 3, 6):
        m = 7*levels
        ell = 1/(32*4**(levels-1))
        s = replacement_horizon(levels)
        failure = 3*m*(2/3)**s
        entropy_error = failure*m*math.log(4) - failure*math.log(failure) - (1-failure)*math.log1p(-failure)
        horizons.append(dict(levels=levels, pairs=m, inverse_smallest_radius=round(1/ell),
                             with_replacement_ticks_per_pair=s, failure_bound=failure,
                             entropy_error_per_smallest_radius_fourth=entropy_error/ell**4,
                             modular_error_per_smallest_radius_fourth=failure*m*math.log(8)/ell**4))
    return dict(schema="oph.one_bit_material_bath.v1", exact_ticks_per_pair=len(SWEEP),
                native_banks=main_banks(), skipped_address=forgotten_address_control(), stochastic_clock=horizons)


def _compare(actual, saved, path=()):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or actual.keys() != saved.keys():
            raise ValueError(f"material bath fields changed at {path}")
        for key, value in actual.items():
            _compare(value, saved[key], path+(key,))
    elif isinstance(actual, list):
        if not isinstance(saved, list) or len(actual) != len(saved):
            raise ValueError(f"material bath record length changed at {path}")
        for i, (value, expected) in enumerate(zip(actual,saved,strict=True)):
            _compare(value, expected, path+(i,))
    elif isinstance(actual, float):
        # Count zeros, probability support and fixed acceptance thresholds are exact.
        numerical_bound = isinstance(path[-1], str) and path[-1].endswith("_less_than")
        if type(saved) not in (int,float) or (actual == 0) != (saved == 0):
            raise ValueError(f"material bath numerical support changed at {path}")
        if numerical_bound:
            matches = actual == saved
        elif "stochastic_clock" in path:
            # The smallest error certificates lie well below ordinary absolute tolerances.
            matches = math.isclose(actual, saved, rel_tol=2e-12, abs_tol=0)
        else:
            matches = math.isclose(actual, saved, rel_tol=2e-12, abs_tol=2e-12)
        if not matches:
            raise ValueError(f"material bath numerical reading changed at {path}")
    elif type(saved) is not type(actual) or actual != saved:
        raise ValueError(f"material bath discrete reading changed at {path}")


def verify(path=None):
    """Replay the instruments and compare every archived record and precision bound."""
    data = deposit()
    saved = HERE/"results.json" if path is None else Path(path)
    _compare(data, json.loads(saved.read_text()))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write:
        (HERE/"results.json").write_text(json.dumps(deposit(), indent=2, allow_nan=False)+"\n")
    else:
        verify()
        print("Replayed every correlated two-memory bath record and the missed-address control.")


if __name__ == "__main__":
    main()
